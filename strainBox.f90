program strain_box
use hdf5
use, intrinsic :: ISO_FORTRAN_ENV
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'
type ragged_array
  !2D array containing vectors of different lengths
  double precision,allocatable::v(:)
  character(len=12) :: iStepChar
end type ragged_array
integer:: ierr,namelen,fileReadErr
logical :: check
integer(kind=int64) :: filesize, expectedSize, dropVol
integer :: im,jm,km
character(len=200) :: namedir, filename
character(len=12) :: iStepChar
character(len=5) :: rankChar
double precision :: modnor
integer,parameter :: maxMom=2,statU=41,posiU=42,veloU=43,MoInU=44,counU=45,topoU=46,voidU=47
integer,parameter :: nxt=256,nyt=256,nzt=256
double precision,parameter :: lx=256,ly=256,lz=256
double precision,parameter :: dx=1,dy=1,dz=1
integer :: rank, ntask
double precision, parameter :: pi=3.14159265358979
!velocity from dataUVW is double precision
integer(kind=int8), dimension(nxt,nyt,nzt) :: top, car
double precision, dimension(nxt,nyt,nzt) :: kur
real, dimension(nxt,nyt,nzt) :: phi
double precision, dimension(nxt,nyt,nzt,3) :: nor, vel
!use int8 for non shared arrays to save memory
integer(kind=int8), dimension(nxt,nyt,nzt) :: topLocal, s_drop, void
integer(kind=int8), dimension(0:1,0:1,0:1) :: paint, neigh
integer :: i,j,k,id,jd,kd,ip,jp,kp,iq,jq,kq,direc,iShifted,mom,last0(3),last1(3),pos(3),ii,jj,error
integer :: cols,paintIt,faceOnCorner,genus,nVoids,nSmallVoids,onInt,handles
double precision, dimension(maxMom,3) :: dropPos
double precision, dimension(maxMom,3) :: dropVel
double precision, dimension(3), parameter :: l = (/lx,ly,lz/)
double precision :: MoI(3,3), work(8), eiVals(3)
double precision :: invSize, sizeMean, invsizeMean, surfConc, surfConcStdDev, diag, deformation, dropArea, dA
double precision :: maxNor, kurMean, kurStdDev, voidSzMn, voidSzStDv, kurInv, kurInvSq
type(ragged_array) :: hist(3) !histogram of drop mass in x, y and z directions
integer(kind=int64) :: nstep, dropCount, dropSize, remain, faces, edges, vertices, voidSize
character(len=200) :: fileEnd 
character(len=100) :: fmtstr
integer, dimension(3) :: decStart, decEnd, nl, decCoord
integer(8), dimension(3), parameter :: decDims = (/1,1,1/), nt = (/nxt,nyt,nzt/)
integer(hid_t) :: file_id, dset_id, filespace, memspace
integer(hsize_t), dimension(3) :: dims, start, count
integer(HSIZE_T), dimension(1) :: dims1d
real(4) :: Cn
complex :: phi_hat(nxt,nyt,nzt)
!type(fftw_plan) :: plan
type(C_PTR)  :: plan
write(*,'(1x,a)') '                 starting number of drops calculation                       '
!create output files
open(42,file='./output/dropCount.dat',form='formatted',position='append')
close(42,status='keep')

! Open the file (read-only)
call h5open_f(error)
call h5fopen_f("/home/alberto.velamartin/drop_time/we_10/run_break_009/field.015.h5", &
                 H5F_ACC_RDONLY_F, file_id, error)
  dims1d=(/1/) 
  dims=(/256,256,256/) 
  call h5dopen_f(file_id, "Cn", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'Cn', Cn
  call h5dopen_f(file_id, "time", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'time', Cn
  call h5dopen_f(file_id, "c", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, phi, dims, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'phi', phi(1,1,1)
  call h5dopen_f(file_id, "We", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
  call h5dclose_f(dset_id, error)
  write(*,*) 'We', Cn
  call h5dclose_f(dset_id, error)
call h5fclose_f(file_id, error)

!plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE)
plan = fftwf_plan_dft_r2c_3d(nxt, nyt, nzt, phi, phi_hat, FFTW_ESTIMATE)

if(error.ne.0) write(6,*) 'fileReadErr'
dropVol = 0 !IC
do k=1,nzt
  do j=1,nyt
    do i=1,nxt
      if(phi(i,j,k).ge.0.5d0)then
        top(i,j,k)=1
        dropVol = dropVol+1 !IC
      else
        top(i,j,k)=0
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
      nor(i,j,k,1)=(phi(ip,j,k)-phi(im,j,k))*(0.5d0/dx)
      nor(i,j,k,2)=(phi(i,jp,k)-phi(i,jm,k))*(0.5d0/dy)
      nor(i,j,k,3)=(phi(i,j,kp)-phi(i,j,km))*(0.5d0/dz)
      !modnor=dsqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2)
      modnor=sqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2)
      !outward pointing normal
      nor(i,j,k,:)=-nor(i,j,k,:)/modnor
    enddo
  enddo
enddo
if(rank.eq.0)then
  open(44,file='./output/dropVol.dat',access='append',form='formatted')
    write(44,'(2i16,E15.7)') nstep,dropVol,dropVol/(nxt*nyt*nzt*1.0)
  close(44,status='keep')
endif
dropVol=0
do k=1,nzt
  do j=1,nyt
    do i=1,nxt
      !find carrier fluid so we can skip it when looking for voids in get_interface. Assume carrier>25%
      if (dropVol.lt.nxt*nyt*nzt/4.and.top(i,j,k).eq.0) then
        car=0
        dropVol=0
        if (rank.eq.0) write(*,*) 'Carrier fraction',1.0*dropVol/nxt/nyt/nzt
      endif
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
      kur(i,j,k)=(nor(ip,j,k,1)-nor(im,j,k,1))*(0.5d0/dx)+ &
                 (nor(i,jp,k,2)-nor(i,jm,k,2))*(0.5d0/dy)+ &
                 (nor(i,j,kp,3)-nor(i,j,km,3))*(0.5d0/dz)
    enddo
  enddo
enddo
filename=trim(namedir)//'dataUVW'//iStepChar
inquire(file=trim(filename),exist=check,size=filesize)
expectedSize = 1.0*nxt*nyt*nzt*4*8 !IC contains u,v,w,p. Cast to real then long integer to prevent overflow
if(rank.eq.0) then
  write(6,*) 'rank ',rank,'reading ',trim(filename)
endif
if((check.eqv..true.).and.(filesize.eq.expectedSize)) then ! input file contains vof double precision
open(67,file=filename,access='stream',status='old',form='unformatted',convert='little_endian')
  read(67) vel
close(67,status='keep')
else
  write(6,*) 'rank ',rank,'wrong file ',trim(filename),'exist ',check,'filesize',filesize,'expect ',expectedSize
  error=1
  return
endif
fileEnd=iStepChar//'rank'//rankChar//'.dat'
open(statU,file='./mpiOutput/statDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(statU,status='keep')
open(posiU,file='./mpiOutput/posiDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(posiU,status='keep')
open(veloU,file='./mpiOutput/veloDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(veloU,status='keep')
open(MoInU,file='./mpiOutput/MoInDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(MoInU,status='keep')
open(topoU,file='./mpiOutput/topoDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(topoU,status='keep')
open(voidU,file='./mpiOutput/voidDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
close(voidU,status='keep')
dropCount=0
sizeMean=0.0
invsizeMean=0.0
surfConc=0.0
surfConcStdDev=0.0
dropPos=0.0
dropVel=0.0
topLocal=top
do direc=1,3
  allocate(hist(direc)%v(nt(direc)))
enddo
hist(1)%iStepChar='1'
hist(2)%iStepChar='2'
hist(3)%iStepChar='3'
!assign a cuboid subset of the domain to the mpi process
decCoord(1) = rank/decDims(2)/decDims(3)
decCoord(2) = rank/decDims(3) - decCoord(1)*decDims(2)
decCoord(3) = rank - decCoord(2)*decDims(3) - decCoord(1)*decDims(3)*decDims(2)
do direc=1,3
  if(mod(nt(direc),decDims(direc)).eq.0) then 
    nl(direc) = nt(direc)/decDims(direc)
    decStart(direc) = decCoord(direc)*nl(direc)+1
    decEnd(direc) = decStart(direc) + nl(direc)-1
  else 
    write(*,*) 'array size not divisible by decDims in direc', direc
    call mpi_finalize(ierr) 
  endif
enddo
do kd=decStart(3), decEnd(3)
  remain=sum(int(topLocal(decStart(1):decEnd(1),&
    decStart(2):decEnd(2),decStart(3):decEnd(3))))
  if(rank.eq.0.and.mod(kd,nl(3)/10).eq.0) write(*,*) 'rank',rank,'start',decStart(3),'kd',kd,'end',decEnd(3),'remain',remain
  do jd=decStart(2), decEnd(2)
    do id=decStart(1), decEnd(1)
      if(topLocal(id,jd,kd).gt.0)then
        !find a single drop with the painter algorithm
        s_drop=0
        s_drop(id,jd,kd)=1
        dropSize=0
        !call flood_fill(id,jd,kd,1,dropSize,s_drop)
        !calculate drop statistics if minimum i j k in the drop lies in the processor's domain
        extremityLoop: do k=1,nzt
          do j=1,nyt
            do i=1,nxt
              if(s_drop(i,j,k).ne.0) exit extremityLoop
            enddo
          enddo
        enddo extremityLoop
        if(i.ge.decStart(1).and.i.le.decEnd(1).and. &
           j.ge.decStart(2).and.j.le.decEnd(2).and. &
           k.ge.decStart(3).and.k.le.decEnd(3)) then
          dropCount = dropCount + 1
          sizeMean = sizeMean + (dropSize - sizeMean) / (1.0*dropCount)
          invSize = 1.0/dropSize
          invsizeMean = invsizeMean + ( invSize - invsizeMean) / (1.0*dropCount)
          dropArea=0.0
          surfConc=0.0
          surfConcStdDev=0.0
          kurMean=0.0d0
          kurStdDev=0.0d0
          kurInv=0.0d0
          kurInvSq=0.0d0
          nVoids=0
          nSmallVoids=0
          voidSzMn=0.0
          voidSzStDv=0.0
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
          dropPos=0.0
          last0=0
          last1=0
          do direc=1,3
            !find beginning and end of drop in x,y,z directions
            do i=1,nt(direc)
              ip=i+1
              if(ip.gt.nt(direc)) ip=ip-nt(direc)
              if(hist(direc)%v(i).gt.0.5*invSize.and.hist(direc)%v(ip).lt.0.5*invSize) last1(direc)=i
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
                !find size of voids
                ip = modulo(i,nxt)+1
                if (s_drop(i,j,k).eq.1.and.s_drop(ip,j,k).eq.0.and.car(ip,j,k).eq.0) then
                  void=0
                  voidSize=0
                  !call flood_fill(ip,j,k,0,voidSize,void)
                  !if (rank.eq.0) 
                  !write(*,*) 'rank',rank,'dropSize',dropSize,'voidSize',voidSize,'sumCar',sum(int(car))
                  do kp=1,nzt
                    do jp=1,nyt
                      do ip=1,nxt
                        if(void(ip,jp,kp).ne.0) car(ip,jp,kp)=1
                      enddo
                    enddo
                  enddo 
                  !count voids smaller than 8 cells separately, because they likely result from vof coalescence
                  if (voidSize.gt.7) then
                    nVoids=nVoids+1
                    voidSzMn = voidSzMn + (voidSize - voidSzMn)/nVoids
                    voidSzStDv = voidSzStDv + (voidSize**2 - voidSzStDv)/nVoids
                  else
                    nSmallVoids=nSmallVoids+1
                  endif
                  open(voidU,file='./mpiOutput/voidDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
                    write(voidU,'(5i16)') rank,dropCount,dropSize,nSmallVoids+nVoids,voidSize
                  close(voidU,status='keep')
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
            write(*,*)'rank',rank,'count',dropCount,'warning: Euler characteristic is odd, v,e,f=',vertices,edges,faces
          endif
          genus = 1-(vertices-edges+faces)/2
          handles=genus+nSmallVoids+nVoids
          kurMean = kurMean / dropArea
          kurStdDev = dsqrt( kurStdDev/dropArea - kurMean**2 )
          kurInv = kurInv / dropArea
          kurInvSq = kurInvSq / dropArea
          surfConc = surfConc/dropArea
          surfConcStdDev = dsqrt(surfConcStdDev/dropArea - surfConc**2)
          voidSzStDv = dsqrt(voidSzStDv - voidSzMn**2)
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
            write(*,*) 'rank',rank,'drop',dropCount,' spans all of domain, deformation is undefined'
            deformation=-1.0
          else
            call DSYEV("N","U",3,MoI,3,eiVals,work,8,error)
            if (error.ne.0) write(*,*) 'dropSize', dropSize, 'eiVals', eiVals
            deformation=sqrt(eiVals(3)/eiVals(1))
          endif
          open(statU,file='./mpiOutput/statDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
            write(statU,'(3i16,6ES16.7E3)') rank,dropCount,dropSize,dropArea,surfConc,surfConcStdDev,deformation,kurMean,kurStdDev
          close(statU,status='keep')
          ! generate format string for writing 
          write(fmtstr,'(a,i0,a)') '(3i16,',maxMom*3,'(ES16.7E3))'
          open(posiU,file='./mpiOutput/posiDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
            write(posiU,fmtstr) rank,dropCount,dropSize,(dropPos(mom,:),mom=1,maxMom)
          close(posiU,status='keep')
          open(veloU,file='./mpiOutput/veloDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
            write(veloU,fmtstr) rank,dropCount,dropSize,(dropVel(mom,:),mom=1,maxMom)
          close(veloU,status='keep')
          open(MoInU,file='./mpiOutput/MoInDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
            write(MoInU,'(3i16,9ES16.7E3,6i16,2ES16.7E3)') rank,dropCount,dropSize,(eiVals(i),i=1,3),&
              MoI(1,1),MoI(1,2),MoI(1,3),MoI(2,2),MoI(2,3),MoI(3,3),last0(:),last1(:),kurInv,kurInvSq
          close(MoInU,status='keep')
          open(topoU,file='./mpiOutput/topoDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
            write(topoU,'(9i16,2ES16.7E3)')rank,dropCount,dropSize,vertices,edges,faces,genus,nSmallVoids,nVoids,voidSzMn,&
              voidSzStDv
          close(topoU,status='keep')
        endif
      endif
    enddo
  enddo
enddo
!2022 Oct 17 Ianto Cannon: inverse of invsizeMean gives the mean drop size with random search
if(dropCount.gt.0) then
  open(counU,file='./output/dropCount.dat',access='append',form='formatted',status='old')
    write(counU,'(4i16,2E15.7)') nstep,rank,dropCount,remain,sizeMean,invsizeMean
  close(counU,status='keep')
endif
write(6,*) 'This is the end'
end program strain_box
